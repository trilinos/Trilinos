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

#ifndef KOKKOS_DEFAULTARITHMETIC_H
#define KOKKOS_DEFAULTARITHMETIC_H

#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Tuple.hpp>
#include <stdexcept>

#include "Kokkos_MultiVector.hpp"
#include "Kokkos_MultiVectorKernelOps.hpp"
#include "Kokkos_NodeHelpers.hpp"
#ifdef HAVE_KOKKOSCLASSIC_TBB
#include "Kokkos_TBBNode.hpp"
#endif
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
#include "Kokkos_OpenMPNode.hpp"
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
#include "Kokkos_TPINode.hpp"
#endif
#ifdef HAVE_KOKKOSCLASSIC_THRUST
#include "Kokkos_ThrustGPUNode.hpp"
#include "cublas.h"
#endif
#include "Kokkos_SerialNode.hpp"
#include <Teuchos_BLAS.hpp>

namespace Kokkos {

  // Class for providing GEMM for a particular Node
  template <typename Scalar, typename Node>
  struct NodeGEMM {
    public:
      static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha, const MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &B, Scalar beta, MultiVector<Scalar,Node> &C) {
        TEUCHOS_TEST_FOR_EXCEPT(true);
      }
  };

  template <typename Scalar>
  struct NodeGEMM<Scalar,SerialNode> {
    public:
      static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha, const MultiVector<Scalar,SerialNode> &A, const MultiVector<Scalar,SerialNode> &B, Scalar beta, MultiVector<Scalar,SerialNode> &C) {
        Teuchos::BLAS<int,Scalar> blas;
        const int m = Teuchos::as<int>(C.getNumRows()),
                  n = Teuchos::as<int>(C.getNumCols()),
                  k = (transA == Teuchos::NO_TRANS ? A.getNumCols() : A.getNumRows()),
                  lda = Teuchos::as<int>(A.getStride()),
                  ldb = Teuchos::as<int>(B.getStride()),
                  ldc = Teuchos::as<int>(C.getStride());
        blas.GEMM(transA, transB, m, n, k, alpha, A.getValues().getRawPtr(), lda, B.getValues().getRawPtr(), ldb, beta, C.getValuesNonConst().getRawPtr(), ldc);
      }
  };

#ifdef HAVE_KOKKOSCLASSIC_TBB
  template <typename Scalar>
  struct NodeGEMM<Scalar,TBBNode> {
    public:
      static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha, const MultiVector<Scalar,TBBNode> &A, const MultiVector<Scalar,TBBNode> &B, Scalar beta, MultiVector<Scalar,TBBNode> &C) {
        Teuchos::BLAS<int,Scalar> blas;
        const int m = Teuchos::as<int>(C.getNumRows()),
                  n = Teuchos::as<int>(C.getNumCols()),
                  k = (transA == Teuchos::NO_TRANS ? A.getNumCols() : A.getNumRows()),
                  lda = Teuchos::as<int>(A.getStride()),
                  ldb = Teuchos::as<int>(B.getStride()),
                  ldc = Teuchos::as<int>(C.getStride());
        blas.GEMM(transA, transB, m, n, k, alpha, A.getValues().getRawPtr(), lda, B.getValues().getRawPtr(), ldb, beta, C.getValuesNonConst().getRawPtr(), ldc);
      }
  };
#endif

#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
  template <typename Scalar>
  struct NodeGEMM<Scalar,TPINode> {
    public:
      static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha, const MultiVector<Scalar,TPINode> &A, const MultiVector<Scalar,TPINode> &B, Scalar beta, MultiVector<Scalar,TPINode> &C) {
#ifndef KOKKOS_DONT_BLOCK_TPI_GEMM
        TPI_Block();
#endif
        Teuchos::BLAS<int,Scalar> blas;
        const int m = Teuchos::as<int>(C.getNumRows()),
                  n = Teuchos::as<int>(C.getNumCols()),
                  k = (transA == Teuchos::NO_TRANS ? A.getNumCols() : A.getNumRows()),
                  lda = Teuchos::as<int>(A.getStride()),
                  ldb = Teuchos::as<int>(B.getStride()),
                  ldc = Teuchos::as<int>(C.getStride());
        blas.GEMM(transA, transB, m, n, k, alpha, A.getValues().getRawPtr(), lda, B.getValues().getRawPtr(), ldb, beta, C.getValuesNonConst().getRawPtr(), ldc);
#ifndef KOKKOS_DONT_BLOCK_TPI_GEMM
        TPI_Unblock();
#endif
      }
  };
#endif

#ifdef HAVE_KOKKOSCLASSIC_OPENMP
  template <typename Scalar>
  struct NodeGEMM<Scalar,OpenMPNode> {
    public:
      static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha, const MultiVector<Scalar,OpenMPNode> &A, const MultiVector<Scalar,OpenMPNode> &B, Scalar beta, MultiVector<Scalar,OpenMPNode> &C) {
        Teuchos::BLAS<int,Scalar> blas;
        const int m = Teuchos::as<int>(C.getNumRows()),
                  n = Teuchos::as<int>(C.getNumCols()),
                  k = (transA == Teuchos::NO_TRANS ? A.getNumCols() : A.getNumRows()),
                  lda = Teuchos::as<int>(A.getStride()),
                  ldb = Teuchos::as<int>(B.getStride()),
                  ldc = Teuchos::as<int>(C.getStride());
        blas.GEMM(transA, transB, m, n, k, alpha, A.getValues().getRawPtr(), lda, B.getValues().getRawPtr(), ldb, beta, C.getValuesNonConst().getRawPtr(), ldc);
      }
  };
#endif

#ifdef HAVE_KOKKOSCLASSIC_THRUST
  template <typename Scalar>
  struct NodeGEMM<Scalar,ThrustGPUNode> {
    public:
      static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha, const MultiVector<Scalar,ThrustGPUNode> &A, const MultiVector<Scalar,ThrustGPUNode> &B, Scalar beta, MultiVector<Scalar,ThrustGPUNode> &C) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "NodeGEMM: ThrustGPUNode has no support for GEMM operations over Scalar=" << Teuchos::typeName(alpha) << ".");
      }
  };


#ifdef HAVE_KOKKOSCLASSIC_CUDA_FLOAT
  template <>
  struct NodeGEMM<float,ThrustGPUNode> {
    public:
      static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, float alpha, const MultiVector<float,ThrustGPUNode> &A, const MultiVector<float,ThrustGPUNode> &B, float beta, MultiVector<float,ThrustGPUNode> &C) {
        const int m = Teuchos::as<int>(C.getNumRows()),
                  n = Teuchos::as<int>(C.getNumCols()),
                  k = (transA == Teuchos::NO_TRANS ? A.getNumCols() : A.getNumRows()),
                  lda = Teuchos::as<int>(A.getStride()),
                  ldb = Teuchos::as<int>(B.getStride()),
                  ldc = Teuchos::as<int>(C.getStride());
        const char char_transA = (transA == Teuchos::NO_TRANS ? 'N' : 'T'),
                   char_transB = (transB == Teuchos::NO_TRANS ? 'N' : 'T');
        cublasSgemm(char_transA, char_transB, m, n, k, alpha, A.getValues().getRawPtr(), lda, B.getValues().getRawPtr(), ldb, beta, C.getValuesNonConst().getRawPtr(), ldc);
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
        cublasStatus info = cublasGetError();
        TEUCHOS_TEST_FOR_EXCEPTION( info != CUBLAS_STATUS_SUCCESS, std::runtime_error, "cublasSgemm failed with status " << info << "." );
#endif
      }
  };
#endif

#ifdef HAVE_KOKKOSCLASSIC_CUDA_DOUBLE
  template <>
  struct NodeGEMM<double,ThrustGPUNode> {
    public:
      static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, double alpha, const MultiVector<double,ThrustGPUNode> &A, const MultiVector<double,ThrustGPUNode> &B, double beta, MultiVector<double,ThrustGPUNode> &C) {
        const int m = Teuchos::as<int>(C.getNumRows()),
                  n = Teuchos::as<int>(C.getNumCols()),
                  k = (transA == Teuchos::NO_TRANS ? A.getNumCols() : A.getNumRows()),
                  lda = Teuchos::as<int>(A.getStride()),
                  ldb = Teuchos::as<int>(B.getStride()),
                  ldc = Teuchos::as<int>(C.getStride());
        const char char_transA = (transA == Teuchos::NO_TRANS ? 'N' : 'T'),
                   char_transB = (transB == Teuchos::NO_TRANS ? 'N' : 'T');
        cublasDgemm(char_transA, char_transB, m, n, k, alpha, A.getValues().getRawPtr(), lda, B.getValues().getRawPtr(), ldb, beta, C.getValuesNonConst().getRawPtr(), ldc);
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
        cublasStatus info = cublasGetError();
        TEUCHOS_TEST_FOR_EXCEPTION( info != CUBLAS_STATUS_SUCCESS, std::runtime_error, "cublasDgemm failed with status " << info << "." );
#endif
      }
  };
#endif

#endif

  /** \brief A traits class providing a generic arithmetic interface for vectors and multivectors.
    */
  template <class MV>
  class DefaultArithmetic {
    // nothing here
  };

  /** \brief Partial specialization of class DefaultArithmetic for MultiVector<Scalar,Node>
    */
  template <class Scalar, class Node>
  class DefaultArithmetic<MultiVector<Scalar,Node> > {

    public:
      //! Initialize multivector to constant value.
      inline static void Init(MultiVector<Scalar,Node> &A, Scalar alpha) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        if (nR*nC == 0) return;
        const size_t stride = A.getStride();
        RCP<Node> node = A.getNode();
        ArrayRCP<Scalar> data = A.getValuesNonConst();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addNonConstBuffer<Scalar>(data);
        rbh.end();
        // prepare op
        InitOp<Scalar> wdp;
        wdp.alpha = alpha;
        if (stride == nR) {
          // one kernel invocation for whole multivector
          wdp.x = data(0,nR*nC).getRawPtr();
          node->template parallel_for<InitOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_t j=0; j<nC; ++j) {
            wdp.x = data(0,nR).getRawPtr();
            node->template parallel_for<InitOp<Scalar> >(0,nR,wdp);
            data += stride;
          }
        }
      }

      //! Set MultiVector to the reciprocal of another: B(i,j) = 1/A(i,j)
      inline static void Recip(MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &B) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        const size_t Bstride = B.getStride();
        TEUCHOS_TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Recip(A,B): A and B must have the same dimensions.");
        RCP<Node> node = B.getNode();
        ArrayRCP<const Scalar> Bdata = B.getValues();
        ArrayRCP<Scalar>       Adata = A.getValuesNonConst();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Bdata);
        rbh.template addNonConstBuffer<Scalar>(Adata);
        rbh.end();
        RecipOp<Scalar> wdp;
        if (A.getStride() == nR && B.getStride() == nR) {
          // one kernel invocation for whole multivector
          wdp.scale = Bdata(0,nR*nC).getRawPtr();
          wdp.x = Adata(0,nR*nC).getRawPtr();
          node->template parallel_for<RecipOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_t j=0; j<nC; ++j) {
            wdp.x = Adata(0,nR).getRawPtr();
            wdp.scale = Bdata(0,nR).getRawPtr();
            node->template parallel_for<RecipOp<Scalar> >(0,nR,wdp);
            Adata += Astride;
            Bdata += Bstride;
          }
        }
      }

      //! Set MultiVector to the scaled element-wise multiple of two others:
      /** C(i,j) = scalarC * C(i,j) + scalarAB * B(i,j) * A(i,1) (A has only 1 column)
       */
      inline static void ElemMult(MultiVector<Scalar,Node> &C, Scalar scalarC,
                               Scalar scalarAB,
                               const MultiVector<Scalar,Node> &A,
                               const MultiVector<Scalar,Node> &B) {
        const size_t nR_A = A.getNumRows();
        const size_t nC_A = A.getNumCols();
        TEUCHOS_TEST_FOR_EXCEPTION(nC_A != 1, std::runtime_error,
                           "DefaultArithmetic<"<<Teuchos::typeName(A) << ">::ElemMult(C,sC,sAB,A,B): A must have just 1 column.");
        const size_t Cstride = C.getStride();
        const size_t Bstride = B.getStride();
        const size_t nC_C = C.getNumCols();
        const size_t nR_C = C.getNumRows();
        TEUCHOS_TEST_FOR_EXCEPTION(nC_C != B.getNumCols() || nR_A != B.getNumRows() || nR_C != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::ElemMult(C,sC,sAB,A,B): A, B and C must have the same num-rows, B and C must have the same num-cols.");
        RCP<Node> node = B.getNode();
        ArrayRCP<Scalar> Cdata = C.getValuesNonConst();
        ArrayRCP<const Scalar> Bdata = B.getValues();
        ArrayRCP<const Scalar>       Adata = A.getValues();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addNonConstBuffer<Scalar>(Cdata);
        rbh.template addConstBuffer<Scalar>(Bdata);
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        MVElemMultOp<Scalar> wdp;
        wdp.scalarX = scalarC;
        wdp.scalarYZ = scalarAB;
        // one kernel invocation for each column
        for (size_t j=0; j<nC_C; ++j) {
          wdp.x = Cdata(0,nR_C).getRawPtr();
          wdp.y = Adata(0,nR_C).getRawPtr();
          wdp.z = Bdata(0,nR_C).getRawPtr();
          node->template parallel_for<MVElemMultOp<Scalar> >(0,nR_C,wdp);
          Cdata += Cstride;
          Bdata += Bstride;
        }
      }

      //! Assign one MultiVector to another
      inline static void Assign(MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &B) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        const size_t Bstride = B.getStride();
        TEUCHOS_TEST_FOR_EXCEPTION(
          nC != B.getNumCols() || nR != B.getNumRows(),
          std::runtime_error,
          "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Assign(A,B): "
          "The MultiVectors A and B do not have the same dimensions.  "
          "A is " << nR << " x " << nC << ", but B is "
          << B.getNumRows() << " x " << B.getNumCols() << ".");
        if (nC*nR == 0) {
          return; // Nothing to do!
        }
        RCP<Node> node = A.getNode();
        ArrayRCP<const Scalar> Bdata = B.getValues();
        ArrayRCP<Scalar>       Adata = A.getValuesNonConst();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Bdata);
        rbh.template addNonConstBuffer<Scalar>(Adata);
        rbh.end();
        // prepare op
        AssignOp<Scalar> wdp;
        if (Astride == nR && Bstride == nR) {
          // one kernel invocation for whole multivector assignment
          wdp.x = Adata(0,nR*nC).getRawPtr();
          wdp.y = Bdata(0,nR*nC).getRawPtr();
          node->template parallel_for<AssignOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_t j=0; j<nC; ++j) {
            wdp.x = Adata(0,nR).getRawPtr();
            wdp.y = Bdata(0,nR).getRawPtr();
            node->template parallel_for<AssignOp<Scalar> >(0,nR,wdp);
            Adata += Astride;
            Bdata += Bstride;
          }
        }
      }

      inline static void Dot(const MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &B, const ArrayView<Scalar> &dots) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        const size_t Bstride = B.getStride();
        TEUCHOS_TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Dot(A,B,dots): A and B must have the same dimensions.");
        TEUCHOS_TEST_FOR_EXCEPTION(nC > Teuchos::as<size_t>(dots.size()), std::runtime_error,
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Dot(A,B,dots): dots must have length as large as number of columns of A and B.");
        if (nR*nC == 0) {
          std::fill( dots.begin(), dots.begin() + nC, Teuchos::ScalarTraits<Scalar>::zero() );
          return;
        }
        RCP<Node> node = A.getNode();
        ArrayRCP<const Scalar> Bdata = B.getValues(),
                                        Adata = A.getValues();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Bdata);
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        DotOp2<Scalar> op;
        for (size_t j=0; j<nC; ++j) {
          op.x = Adata(0,nR).getRawPtr();
          op.y = Bdata(0,nR).getRawPtr();
          dots[j] = node->parallel_reduce(0,nR,op);
          Adata += Astride;
          Bdata += Bstride;
        }
      }

      inline static Scalar Dot(const MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &B) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        TEUCHOS_TEST_FOR_EXCEPTION(nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Dot(A,B,dots): A and B must have the same number of rows.");
        if (nR*nC == 0) {
          return Teuchos::ScalarTraits<Scalar>::zero();
        }
        RCP<Node> node = A.getNode();
        ArrayRCP<const Scalar> Bdata = B.getValues(0),
                                        Adata = A.getValues(0);
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Bdata);
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        DotOp2<Scalar> op;
        op.x = Adata(0,nR).getRawPtr();
        op.y = Bdata(0,nR).getRawPtr();
        return node->parallel_reduce(0,nR,op);
      }

      inline static void GEMM(MultiVector<Scalar,Node> &C, Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha, const MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &B, Scalar beta) {
        NodeGEMM<Scalar,Node>::GEMM(transA, transB, alpha, A, B, beta, C);
      }

      inline static void GESUM(MultiVector<Scalar,Node> &B, Scalar alpha, const MultiVector<Scalar,Node> &A, Scalar beta) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        const size_t Bstride = B.getStride();
        TEUCHOS_TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::GESUM(B,alpha,A,beta): A and B must have the same dimensions.");
        RCP<Node> node = B.getNode();
        ArrayRCP<const Scalar> Adata = A.getValues();
        ArrayRCP<Scalar>       Bdata = B.getValuesNonConst();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.template addNonConstBuffer<Scalar>(Bdata);
        rbh.end();
        GESUMOp<Scalar> wdp;
        wdp.alpha = alpha;
        wdp.beta  = beta;
        if (Astride == nR && Bstride == nR) {
          // one kernel invocation for whole multivector
          wdp.y = Bdata(0,nR*nC).getRawPtr();
          wdp.x = Adata(0,nR*nC).getRawPtr();
          node->template parallel_for<GESUMOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_t j=0; j<nC; ++j) {
            wdp.y = Bdata(0,nR).getRawPtr();
            wdp.x = Adata(0,nR).getRawPtr();
            node->template parallel_for<GESUMOp<Scalar> >(0,nR,wdp);
            Adata += Astride;
            Bdata += Bstride;
          }
        }
      }

      inline static void GESUM(MultiVector<Scalar,Node> &C, Scalar alpha, const MultiVector<Scalar,Node> &A, Scalar beta, const MultiVector<Scalar,Node> &B, Scalar gamma) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        const size_t Bstride = B.getStride();
        const size_t Cstride = C.getStride();
        TEUCHOS_TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::GESUM(C,alpha,A,beta,B,gamma): A and B must have the same dimensions.");
        RCP<Node> node = B.getNode();
        ArrayRCP<const Scalar> Adata = A.getValues(),
                                        Bdata = B.getValues();
        ArrayRCP<Scalar>       Cdata = C.getValuesNonConst();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.template addConstBuffer<Scalar>(Bdata);
        rbh.template addNonConstBuffer<Scalar>(Cdata);
        rbh.end();
        GESUMOp3<Scalar> wdp;
        wdp.alpha = alpha;
        wdp.beta  = beta;
        wdp.gamma = gamma;
        if (Astride == nR && Bstride == nR && Cstride == nR) {
          // one kernel invocation for whole multivector
          wdp.z = Cdata(0,nR*nC).getRawPtr();
          wdp.y = Bdata(0,nR*nC).getRawPtr();
          wdp.x = Adata(0,nR*nC).getRawPtr();
          node->template parallel_for<GESUMOp3<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_t j=0; j<nC; ++j) {
            wdp.z = Cdata(0,nR).getRawPtr();
            wdp.y = Bdata(0,nR).getRawPtr();
            wdp.x = Adata(0,nR).getRawPtr();
            node->template parallel_for<GESUMOp3<Scalar> >(0,nR,wdp);
            Adata += Astride;
            Bdata += Bstride;
            Cdata += Cstride;
          }
        }
      }

      inline static void Norm1(const MultiVector<Scalar,Node> &A, const ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        TEUCHOS_TEST_FOR_EXCEPTION(nC > Teuchos::as<size_t>(norms.size()), std::runtime_error,
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Norm1(A,norms): norms must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          std::fill( norms.begin(), norms.begin() + nC, Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero() );
          return;
        }
        RCP<Node> node = A.getNode();
        ArrayRCP<const Scalar> Adata = A.getValues();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        SumAbsOp<Scalar> op;
        for (size_t j=0; j<nC; ++j) {
          op.x = Adata(0,nR).getRawPtr();
          norms[j] = node->parallel_reduce(0,nR,op);
          Adata += Astride;
        }
      }

      inline static typename Teuchos::ScalarTraits<Scalar>::magnitudeType
      Norm1(const MultiVector<Scalar,Node> &A) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        if (nR*nC == 0) {
          return Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();
        }
        RCP<Node> node = A.getNode();
        ArrayRCP<const Scalar> Adata = A.getValues(0);
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        SumAbsOp<Scalar> op;
        op.x = Adata(0,nR).getRawPtr();
        return node->parallel_reduce(0,nR,op);
      }

      inline static void Sum(const MultiVector<Scalar,Node> &A, const ArrayView<Scalar> &sums) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        TEUCHOS_TEST_FOR_EXCEPTION(nC > (size_t)sums.size(), std::runtime_error,
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Sum(A,sums): sums must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          std::fill( sums.begin(), sums.begin() + nC, Teuchos::ScalarTraits<Scalar>::zero() );
          return;
        }
        RCP<Node> node = A.getNode();
        ArrayRCP<const Scalar> Adata = A.getValues();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        SumOp<Scalar> op;
        for (size_t j=0; j<nC; ++j) {
          op.x = Adata(0,nR).getRawPtr();
          sums[j] = node->parallel_reduce(0,nR,op);
          Adata += Astride;
        }
      }

      inline static Scalar Sum(const MultiVector<Scalar,Node> &A) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        if (nR*nC == 0) {
          return Teuchos::ScalarTraits<Scalar>::zero();
        }
        RCP<Node> node = A.getNode();
        ArrayRCP<const Scalar> Adata = A.getValues(0);
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        SumOp<Scalar> op;
        op.x = Adata(0,nR).getRawPtr();
        return node->parallel_reduce(0,nR,op);
      }

      inline static typename Teuchos::ScalarTraits<Scalar>::magnitudeType NormInf(const MultiVector<Scalar,Node> &A) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        if (nR*nC == 0) {
          return Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();
        }
        RCP<Node> node = A.getNode();
        ArrayRCP<const Scalar> Adata = A.getValues(0);
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        MaxAbsOp<Scalar> op;
        op.x = Adata(0,nR).getRawPtr();
        return node->parallel_reduce(0,nR,op);
      }

      inline static void NormInf(const MultiVector<Scalar,Node> &A, const ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        TEUCHOS_TEST_FOR_EXCEPTION(nC > Teuchos::as<size_t>(norms.size()), std::runtime_error,
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::NormInf(A,norms): norms must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          std::fill( norms.begin(), norms.begin() + nC, Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero() );
          return;
        }
        RCP<Node> node = A.getNode();
        ArrayRCP<const Scalar> Adata = A.getValues();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        MaxAbsOp<Scalar> op;
        for (size_t j=0; j<nC; ++j) {
          op.x = Adata(0,nR).getRawPtr();
          norms[j] = node->parallel_reduce(0,nR,op);
          Adata += Astride;
        }
      }

      inline static void Norm2Squared(const MultiVector<Scalar,Node> &A, const ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        TEUCHOS_TEST_FOR_EXCEPTION(nC > Teuchos::as<size_t>(norms.size()), std::runtime_error,
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Norm2Squared(A,norms): norms must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          std::fill( norms.begin(), norms.begin() + nC, Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero() );
          return;
        }
        RCP<Node> node = A.getNode();
        ArrayRCP<const Scalar> Adata = A.getValues();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        DotOp1<Scalar> op;
        for (size_t j=0; j<nC; ++j) {
          op.x = Adata(0,nR).getRawPtr();
          norms[j] = node->parallel_reduce(0,nR,op);
          Adata += Astride;
        }
      }

      inline static typename Teuchos::ScalarTraits<Scalar>::magnitudeType
      Norm2Squared(const MultiVector<Scalar,Node> &A) {
        const size_t nR = A.getNumRows();
        if (nR == 0) {
          return Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();
        }
        RCP<Node> node = A.getNode();
        ArrayRCP<const Scalar> Adata = A.getValues(0);
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        DotOp1<Scalar> op;
        op.x = Adata(0,nR).getRawPtr();
        return node->parallel_reduce(0,nR,op);
      }

      inline static typename Teuchos::ScalarTraits<Scalar>::magnitudeType
      WeightedNorm(const MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &weightVector) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        if (nR*nC == 0) {
          return Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero();
        }
        RCP<Node> node = A.getNode();
        ArrayRCP<const Scalar> Adata = A.getValues(0),
                                        Wdata = weightVector.getValues(0);
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.template addConstBuffer<Scalar>(Wdata);
        rbh.end();
        WeightNormOp<Scalar> op;
        op.x = Adata(0,nR).getRawPtr();
        op.w = Wdata(0,nR).getRawPtr();
        return node->parallel_reduce(0,nR,op);
      }

      inline static void WeightedNorm(const MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &weightVector, const ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride(),
                     Wstride = weightVector.getStride();
        TEUCHOS_TEST_FOR_EXCEPTION(nC > Teuchos::as<size_t>(norms.size()), std::runtime_error,
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Norm1(A,norms): norms must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          std::fill( norms.begin(), norms.begin() + nC, Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero() );
          return;
        }
        RCP<Node> node = A.getNode();
        ArrayRCP<const Scalar> Adata = A.getValues(),
                                        Wdata = weightVector.getValues();
        const bool OneW = (weightVector.getNumCols() == 1);
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.template addConstBuffer<Scalar>(Wdata);
        rbh.end();
        WeightNormOp<Scalar> op;
        if (OneW) {
          op.w = Wdata(0,nR).getRawPtr();
          for (size_t j=0; j<nC; ++j) {
            op.x = Adata(0,nR).getRawPtr();
            norms[j] = node->parallel_reduce(0,nR,op);
            Adata += Astride;
          }
        }
        else {
          for (size_t j=0; j<nC; ++j) {
            op.x = Adata(0,nR).getRawPtr();
            op.w = Wdata(0,nR).getRawPtr();
            norms[j] = node->parallel_reduce(0,nR,op);
            Adata += Astride;
            Wdata += Wstride;
          }
        }
      }

      inline static void Random(MultiVector<Scalar,Node> &A) {
        // TODO: consider adding rand() functionality to node
        // in the meantime, just generate random numbers via Teuchos and then copy to node
        typedef Teuchos::ScalarTraits<Scalar> SCT;
        const size_t stride = A.getStride();
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        if (nR*nC == 0) return;
        RCP<Node> node = A.getNode();
        ArrayRCP<Scalar> Adata = A.getValuesNonConst();
        // we'll overwrite all data covered by the multivector, but not off-stride data
        // therefore, we are write-only only in the case that stride=nR
        ReadWriteOption rw = (stride == nR ? WriteOnly : ReadWrite);
        ArrayRCP<Scalar> mvdata = node->template viewBufferNonConst<Scalar>(rw,stride*(nC-1)+nR,Adata);
        for (size_t j=0; j<nC; ++j) {
          for (size_t i=0; i<nR; ++i) {
            mvdata[j*stride + i] = SCT::random();
          }
        }
        mvdata = null;
      }

      inline static void Abs(MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &B) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        const size_t Bstride = B.getStride();
        TEUCHOS_TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Abs(A,B): A and B must have the same dimensions.");
        if (nC*nR == 0) return;
        RCP<Node> node = A.getNode();
        ArrayRCP<const Scalar> Bdata = B.getValues();
        ArrayRCP<Scalar>       Adata = A.getValuesNonConst();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Bdata);
        rbh.template addNonConstBuffer<Scalar>(Adata);
        rbh.end();
        // prepare op
        AbsOp<Scalar> wdp;
        if (Astride == nR && Bstride == nR) {
          // one kernel invocation for whole multivector assignment
          wdp.x = Adata(0,nR*nC).getRawPtr();
          wdp.y = Bdata(0,nR*nC).getRawPtr();
          node->template parallel_for<AbsOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_t j=0; j<nC; ++j) {
            wdp.x = Adata(0,nR).getRawPtr();
            wdp.y = Bdata(0,nR).getRawPtr();
            node->template parallel_for<AbsOp<Scalar> >(0,nR,wdp);
            Adata += Astride;
            Bdata += Bstride;
          }
        }
      }

      inline static void Scale(MultiVector<Scalar,Node> &B, Scalar alpha, const MultiVector<Scalar,Node> &A) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        const size_t Bstride = B.getStride();
        TEUCHOS_TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Scale(B,alpha,A): A and B must have the same dimensions.");
        RCP<Node> node = B.getNode();
        ArrayRCP<const Scalar> Adata = A.getValues();
        ArrayRCP<Scalar>       Bdata = B.getValuesNonConst();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.template addNonConstBuffer<Scalar>(Bdata);
        rbh.end();
        MVScaleOp<Scalar> wdp;
        wdp.alpha = alpha;
        if (Astride == nR && Bstride == nR) {
          // one kernel invocation for whole multivector
          wdp.x = Bdata(0,nR*nC).getRawPtr();
          wdp.y = Adata(0,nR*nC).getRawPtr();
          node->template parallel_for<MVScaleOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_t j=0; j<nC; ++j) {
            wdp.x = Bdata(0,nR).getRawPtr();
            wdp.y = Adata(0,nR).getRawPtr();
            node->template parallel_for<MVScaleOp<Scalar> >(0,nR,wdp);
            Adata += Astride;
            Bdata += Bstride;
          }
        }
      }

      inline static void Scale(MultiVector<Scalar,Node> &A, Scalar alpha) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t stride = A.getStride();
        RCP<Node> node = A.getNode();
        ArrayRCP<Scalar> data = A.getValuesNonConst();
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addNonConstBuffer<Scalar>(data);
        rbh.end();
        // prepare op
        SingleScaleOp<Scalar> wdp;
        wdp.alpha = alpha;
        if (stride == nR) {
          // one kernel invocation for whole multivector
          wdp.x = data(0,nR*nC).getRawPtr();
          node->template parallel_for<SingleScaleOp<Scalar> >(0,nR*nC,wdp);
        }
        else {
          // one kernel invocation for each column
          for (size_t j=0; j<nC; ++j) {
            wdp.x = data(0,nR).getRawPtr();
            node->template parallel_for<SingleScaleOp<Scalar> >(0,nR,wdp);
            data += stride;
          }
        }
      }

      inline static void initializeValues(MultiVector<Scalar,Node> &A,
                                   size_t numRows, size_t numCols,
                                   const ArrayRCP<Scalar> &values,
                                   size_t stride) {
        A.initializeValues(numRows,numCols,values,stride);
      }

      inline static ArrayRCP<const Scalar> getValues(const MultiVector<Scalar,Node> &A) {
        return A.getValues();
      }

      inline static ArrayRCP<const Scalar> getValues(const MultiVector<Scalar,Node> &A, size_t j) {
        return A.getValues(j);
      }

      inline static ArrayRCP<Scalar> getValuesNonConst(MultiVector<Scalar,Node> &A) {
        return A.getValuesNonConst();
      }

      inline static ArrayRCP<Scalar> getValuesNonConst(MultiVector<Scalar,Node> &A, size_t j) {
        return A.getValuesNonConst(j);
      }

      inline static size_t getNumRows(const MultiVector<Scalar,Node> &A) {
        return A.getNumRows();
      }

      inline static size_t getNumCols(const MultiVector<Scalar,Node> &A) {
        return A.getNumCols();
      }

      inline static size_t getStride(const MultiVector<Scalar,Node> &A) {
        return A.getStride();
      }

      inline static RCP<Node> getNode(const MultiVector<Scalar,Node> &A) {
        return A.getNode();
      }

  };

}

#endif
