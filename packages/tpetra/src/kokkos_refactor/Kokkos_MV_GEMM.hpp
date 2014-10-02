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

//Note this code lives only temporarily in Tpetra
//As soon as GEMM kernels exist in KokkosLinAlg and thus a depnedency on Teuchos
//can be eliminated the code will move to KokkosLinAlg.

#if defined(KOKKOS_MULTIVECTOR_H_) && defined(TPETRA_KOKKOS_REFACTOR_MULTIVECTOR_DEF_HPP)

#include<Teuchos_BLAS.hpp>
#include<Teuchos_as.hpp>

#ifdef KOKKOS_HAVE_CUDA
#include<cublas.h>
#endif
namespace Kokkos {
  namespace Impl {

    template<class ViewType>
    size_t getStride2DView (ViewType A) {
      size_t stride[8];
      A.stride (stride);
      return A.dimension_1 () > 1 ? stride[1] : A.dimension_0 ();
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
          const Scalar alpha,
          View<const Scalar**, LayoutLeft, DeviceType> A,
          View<const Scalar**, LayoutLeft, DeviceType> B,
          const Scalar beta,
          View<Scalar**, LayoutLeft, DeviceType> C)
    {
      Teuchos::BLAS<int,Scalar> blas;
      const int m = static_cast<int> (C.dimension_0 ()),
        n = static_cast<int> (C.dimension_1 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_1 () : A.dimension_0 ()),
        lda = static_cast<int> (Impl::getStride2DView (A)),
        ldb = static_cast<int> (Impl::getStride2DView (B)),
        ldc = static_cast<int> (Impl::getStride2DView (C));
      // For some BLAS implementations (e.g., MKL), GEMM when B has
      // one column may be signficantly less efficient than GEMV.
      if (n == 1 && transB == Teuchos::NO_TRANS) {
        blas.GEMV (transA, A.dimension_0 (), A.dimension_1 (), alpha,
                   A.ptr_on_device(), lda,
                   B.ptr_on_device(), static_cast<int> (1),
                   beta, C.ptr_on_device(), static_cast<int> (1));
      }
      else {
        blas.GEMM (transA, transB, m, n, k, alpha,
                   A.ptr_on_device(), lda,
                   B.ptr_on_device(), ldb,
                   beta, C.ptr_on_device(), ldc);
      }
    }
  };

//   template <typename Scalar>
//   struct DeviceGEMM<Scalar,Serial> {
//     public:
//       static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha,
//           View<const Scalar**,LayoutLeft,Serial> A, View<const Scalar**,LayoutLeft,Serial> B,
//           Scalar beta, View<Scalar**,Serial> C) {
//         Teuchos::BLAS<int,Scalar> blas;
//         const int m = Teuchos::as<int>(C.dimension_0()),
//                   n = Teuchos::as<int>(C.dimension_1()),
//                   k = (transA == Teuchos::NO_TRANS ? A.dimension_1() : A.dimension_0()),
//                   lda = Teuchos::as<int>(Impl::getStride2DView(A)),
//                   ldb = Teuchos::as<int>(Impl::getStride2DView(B)),
//                   ldc = Teuchos::as<int>(Impl::getStride2DView(C));
//         // For some BLAS implementations (i.e. MKL), GEMM when B has one column
//         // is signficantly less efficient
//         if (n == 1 && transB == Teuchos::NO_TRANS)
//           blas.GEMV(transA, A.dimension_0(), A.dimension_1(), alpha, A.ptr_on_device(), lda, B.ptr_on_device(), Teuchos::as<int>(1), beta, C.ptr_on_device(), Teuchos::as<int>(1));
//         else
//           blas.GEMM(transA, transB, m, n, k, alpha, A.ptr_on_device(), lda, B.ptr_on_device(), ldb, beta, C.ptr_on_device(), ldc);
//       }
//   };

// #ifdef KOKKOS_HAVE_PTHREAD
//   template <typename Scalar>
//   struct DeviceGEMM<Scalar,Threads> {
//     public:
//       static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha,
//           View<const Scalar**,LayoutLeft,Threads> A, View<const Scalar**,LayoutLeft,Threads> B,
//           Scalar beta, View<Scalar**,LayoutLeft,Threads> C) {
//         Teuchos::BLAS<int,Scalar> blas;
//         const int m = Teuchos::as<int>(C.dimension_0()),
//                   n = Teuchos::as<int>(C.dimension_1()),
//                   k = (transA == Teuchos::NO_TRANS ? A.dimension_1() : A.dimension_0()),
//                   lda = Teuchos::as<int>(Impl::getStride2DView(A)),
//                   ldb = Teuchos::as<int>(Impl::getStride2DView(B)),
//                   ldc = Teuchos::as<int>(Impl::getStride2DView(C));
//         blas.GEMM(transA, transB, m, n, k, alpha, A.ptr_on_device(), lda, B.ptr_on_device(), ldb, beta, C.ptr_on_device(), ldc);
//       }
//   };
// #endif

// #ifdef KOKKOS_HAVE_OPENMP
//   template <typename Scalar>
//   struct DeviceGEMM<Scalar,OpenMP> {
//     public:
//       static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha,
//           View<const Scalar**,LayoutLeft,OpenMP> A, View<const Scalar**,LayoutLeft,OpenMP> B,
//           Scalar beta, View<Scalar**,LayoutLeft,OpenMP> C) {
//         Teuchos::BLAS<int,Scalar> blas;
//         const int m = Teuchos::as<int>(C.dimension_0()),
//                   n = Teuchos::as<int>(C.dimension_1()),
//                   k = (transA == Teuchos::NO_TRANS ? A.dimension_1() : A.dimension_0()),
//                   lda = Teuchos::as<int>(Impl::getStride2DView(A)),
//                   ldb = Teuchos::as<int>(Impl::getStride2DView(B)),
//                   ldc = Teuchos::as<int>(Impl::getStride2DView(C));
//         blas.GEMM(transA, transB, m, n, k, alpha, A.ptr_on_device(), lda, B.ptr_on_device(), ldb, beta, C.ptr_on_device(), ldc);
//       }
//   };
// #endif

#ifdef KOKKOS_HAVE_CUDA
  template <typename Scalar>
  struct DeviceGEMM<Scalar,Cuda> {
    public:
      static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha,
          View<const Scalar**,LayoutLeft,Cuda> A, View<const Scalar**,LayoutLeft,Cuda> B,
          Scalar beta, View<Scalar**,LayoutLeft,Cuda> C) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "DeviceGEMM: Kokkos::Cuda has no support for GEMM operations over Scalar=" << Teuchos::typeName(alpha) << ".");
      }
  };


  template <>
  struct DeviceGEMM<float,Cuda> {
    public:
      static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, float alpha,
          View<const float**,LayoutLeft,Cuda> A, View<const float**,LayoutLeft,Cuda> B,
          float beta, View<float**,LayoutLeft,Cuda> C) {
        const int m = Teuchos::as<int>(C.dimension_0()),
                  n = Teuchos::as<int>(C.dimension_1()),
                  k = (transA == Teuchos::NO_TRANS ? A.dimension_1() : A.dimension_0()),
                  lda = Teuchos::as<int>(Impl::getStride2DView(A)),
                  ldb = Teuchos::as<int>(Impl::getStride2DView(B)),
                  ldc = Teuchos::as<int>(Impl::getStride2DView(C));
        const char char_transA = (transA == Teuchos::NO_TRANS ? 'N' : 'T'),
                   char_transB = (transB == Teuchos::NO_TRANS ? 'N' : 'T');
        cublasSgemm(char_transA, char_transB, m, n, k, alpha, A.ptr_on_device(), lda, B.ptr_on_device(), ldb, beta, C.ptr_on_device(), ldc);
#ifdef HAVE_KOKKOS_DEBUG
        cublasStatus info = cublasGetError();
        TEUCHOS_TEST_FOR_EXCEPTION( info != CUBLAS_STATUS_SUCCESS, std::runtime_error, "cublasSgemm failed with status " << info << "." );
#endif
      }
  };

  template <>
  struct DeviceGEMM<double,Cuda> {
    public:
      static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, double alpha,
          View<const double**,LayoutLeft,Cuda> A, View<const double**,LayoutLeft,Cuda> B,
          double beta, View<double**,LayoutLeft,Cuda> C) {
        const int m = Teuchos::as<int>(C.dimension_0()),
                  n = Teuchos::as<int>(C.dimension_1()),
                  k = (transA == Teuchos::NO_TRANS ? A.dimension_1() : A.dimension_0()),
                  lda = Teuchos::as<int>(Impl::getStride2DView(A)),
                  ldb = Teuchos::as<int>(Impl::getStride2DView(B)),
                  ldc = Teuchos::as<int>(Impl::getStride2DView(C));
        const char char_transA = (transA == Teuchos::NO_TRANS ? 'N' : 'T'),
                   char_transB = (transB == Teuchos::NO_TRANS ? 'N' : 'T');
        cublasDgemm(char_transA, char_transB, m, n, k, alpha, A.ptr_on_device(), lda, B.ptr_on_device(), ldb, beta, C.ptr_on_device(), ldc);
#ifdef HAVE_KOKKOS_DEBUG
        cublasStatus info = cublasGetError();
        TEUCHOS_TEST_FOR_EXCEPTION( info != CUBLAS_STATUS_SUCCESS, std::runtime_error, "cublasDgemm failed with status " << info << "." );
#endif
      }
  };


#endif
}
#endif

