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

#ifndef KOKKOS_DEFAULTARITHMETIC_H
#define KOKKOS_DEFAULTARITHMETIC_H

#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Tuple.hpp>
#include <stdexcept>

#include "Kokkos_MultiVector.hpp"
#include "Kokkos_NodeHelpers.hpp"
#ifdef HAVE_KOKKOS_TBB
#include "Kokkos_TBBNode.hpp"
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
#include "Kokkos_TPINode.hpp"
#endif
#include "Kokkos_SerialNode.hpp"
#include <Teuchos_BLAS.hpp>

#ifndef KERNEL_PREFIX 
  #define KERNEL_PREFIX
#endif

namespace Kokkos {

  template <class Scalar>
  struct InitOp {
    Scalar *x;
    Scalar alpha;
    inline KERNEL_PREFIX void execute(size_t i) const
    {
      x[i] = alpha;
    }
  };

  template <class Scalar>
  struct AssignOp {
    Scalar *x;
    const Scalar *y;
    inline KERNEL_PREFIX void execute(size_t i) const
    {
      x[i] = y[i];
    }
  };

  template <class Scalar>
  struct SingleScaleOp {
    Scalar alpha;
    Scalar *x;
    inline KERNEL_PREFIX void execute(size_t i) const
    {
      Scalar tmp = x[i];
      x[i] = alpha*tmp;
    }
  };

  template <class Scalar>
  struct MVScaleOp {
    Scalar alpha;
    const Scalar *y;
    Scalar *x;
    inline KERNEL_PREFIX void execute(size_t i) const
    {
      x[i] = alpha*y[i];
    }
  };

  template <class Scalar>
  struct AbsOp {
    const Scalar *y;
    Scalar *x;
    inline KERNEL_PREFIX void execute(size_t i) const
    {
      x[i] = Teuchos::ScalarTraits<Scalar>::magnitude(y[i]);
    }
  };

  template <class Scalar>
  struct RecipOp {
    const Scalar *scale;
    Scalar *x;
    inline KERNEL_PREFIX void execute(size_t i) const
    {
      x[i] = Teuchos::ScalarTraits<Scalar>::one() / scale[i];
    }
  };

  template <class Scalar>
  struct GESUMOp {
    const Scalar *x;
    Scalar *y;
    Scalar alpha, beta;
    inline KERNEL_PREFIX void execute(size_t i) const
    {
      Scalar tmp = y[i];
      y[i] = alpha * x[i] + beta * tmp;
    }
  };

  template <class Scalar>
  struct GESUMOp3 {
    const Scalar *x, *y;
    Scalar *z;
    Scalar alpha, beta, gamma;
    inline KERNEL_PREFIX void execute(size_t i) const
    {
      Scalar tmp = z[i];
      z[i] = alpha * x[i] + beta * y[i] + gamma * tmp;
    }
  };

  template <class Scalar>
  struct SumAbsOp {
    typedef  Teuchos::ScalarTraits<Scalar> SCT;
    typedef  typename SCT::magnitudeType   Magnitude;
    const Scalar *x;
    typedef  Magnitude ReductionType;
    inline static Magnitude identity() {return Teuchos::ScalarTraits<Magnitude>::zero();}
    Magnitude reduce(Magnitude x, Magnitude y) {return x+y;}
    Magnitude generate(size_t i) {
      return SCT::magnitude(x[i]);
    }
  };

  template <class Scalar>
  struct WeightNormOp {
    typedef  Teuchos::ScalarTraits<Scalar> SCT;
    typedef  typename SCT::magnitudeType   Magnitude;
    const Scalar *x, *w;
    typedef  Magnitude ReductionType;
    inline static Magnitude identity() {return Teuchos::ScalarTraits<Magnitude>::zero();}
    Magnitude reduce(Magnitude x, Magnitude y) {return x+y;}
    Magnitude generate(size_t i) {
      Scalar tmp = x[i] / w[i];
      return SCT::real( SCT::conjugate(tmp)*tmp );
    }
  };

  template <class Scalar>
  struct SumOp {
    const Scalar *x;
    typedef  Scalar ReductionType;
    inline static Scalar identity() {return Teuchos::ScalarTraits<Scalar>::zero();}
    Scalar reduce(Scalar x, Scalar y) {return x+y;}
    Scalar generate(size_t i) { return x[i]; }
  };

  template <class Scalar>
  struct MaxAbsOp {
    typedef  Teuchos::ScalarTraits<Scalar> SCT;
    typedef  typename SCT::magnitudeType   Magnitude;
    const Scalar *x;
    typedef  Magnitude ReductionType;
    inline static Magnitude identity() {return Teuchos::ScalarTraits<Magnitude>::zero();}
    Magnitude reduce(Magnitude x, Magnitude y) {return std::max(x,y);}
    Magnitude generate(size_t i) {
      return SCT::magnitude(x[i]);
    }
  };

  template <class Scalar>
  struct DotOp1 {
    typedef  Teuchos::ScalarTraits<Scalar> SCT;
    typedef  typename SCT::magnitudeType   Magnitude;
    const Scalar *x;
    typedef  Magnitude ReductionType;
    inline static Magnitude identity() {return Teuchos::ScalarTraits<Magnitude>::zero();}
    Magnitude reduce(Magnitude x, Magnitude y) {return x+y;}
    Magnitude generate(size_t i) {
      Scalar xi = x[i]; 
      return SCT::real( SCT::conjugate(xi)*xi );
    }
  };

  template <class Scalar>
  struct DotOp2 {
    typedef Teuchos::ScalarTraits<Scalar> SCT;
    const Scalar *x, *y;
    typedef  Scalar ReductionType;
    inline static Scalar identity() {return SCT::zero();}
    Scalar reduce(Scalar x, Scalar y) {return x+y;}
    Scalar generate(size_t i) {
      Scalar xi = x[i]; Scalar yi = y[i];
      return SCT::real( SCT::conjugate(xi)*yi );
    }
  };

  //! Class for providing GEMM for a particular Node
  template <class Scalar, class Node> 
  class NodeGEMM {
    public:
      static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha, const MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &B, Scalar beta, MultiVector<Scalar,Node> &C) {
        TEST_FOR_EXCEPT(true);
      }
  };

  template <class Scalar> 
  class NodeGEMM<Scalar,SerialNode> {
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

#ifdef HAVE_KOKKOS_TBB
  template <class Scalar> 
  class NodeGEMM<Scalar,TBBNode> {
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

#ifdef HAVE_KOKKOS_THREADPOOL
  template <class Scalar> 
  class NodeGEMM<Scalar,TPINode> {
    public:
      static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha, const MultiVector<Scalar,TPINode> &A, const MultiVector<Scalar,TPINode> &B, Scalar beta, MultiVector<Scalar,TPINode> &C) {
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

  //! Class DefaultArithmetic, unimplemented
  template <class MV>
  class DefaultArithmetic {
    // nothing here
  };

  //! partial specialization of class DefaultArithmetic, for Kokkos::MultiVector<Scalar,Node>
  template <class Scalar, class Node>
  class DefaultArithmetic<MultiVector<Scalar,Node> > {

    public:
      //! Initialize multivector to constant value.
      inline static void Init(MultiVector<Scalar,Node> &A, Scalar alpha) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t stride = A.getStride();
        Teuchos::RCP<Node> node = A.getNode();
        Teuchos::ArrayRCP<Scalar> data = A.getValuesNonConst();
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
        TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Recip(A,B): A and B must have the same dimensions.");
        Teuchos::RCP<Node> node = B.getNode();
        Teuchos::ArrayRCP<const Scalar> Bdata = B.getValues();
        Teuchos::ArrayRCP<Scalar>       Adata = A.getValuesNonConst();
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

      //! Assign one MultiVector to another
      inline static void Assign(MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &B) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        const size_t Bstride = B.getStride();
        TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Assign(A,B): A and B must have the same dimensions.");
        if (nC*nR == 0) return;
        Teuchos::RCP<Node> node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Bdata = B.getValues();
        Teuchos::ArrayRCP<Scalar>       Adata = A.getValuesNonConst();
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

      inline static void Dot(const MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &B, const Teuchos::ArrayView<Scalar> &dots) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        const size_t Bstride = B.getStride();
        TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Dot(A,B,dots): A and B must have the same dimensions.");
        TEST_FOR_EXCEPTION(nC > Teuchos::as<size_t>(dots.size()), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Dot(A,B,dots): dots must have length as large as number of columns of A and B.");
        if (nR*nC == 0) {
          std::fill( dots.begin(), dots.begin() + nC, Teuchos::ScalarTraits<Scalar>::zero() );
          return;
        }
        Teuchos::RCP<Node> node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Bdata = B.getValues(),
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
        TEST_FOR_EXCEPTION(nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Dot(A,B,dots): A and B must have the same number of rows.");
        if (nR*nC == 0) {
          return Teuchos::ScalarTraits<Scalar>::zero();
        }
        Teuchos::RCP<Node> node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Bdata = B.getValues(0),
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
        TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::GESUM(B,alpha,A,beta): A and B must have the same dimensions.");
        Teuchos::RCP<Node> node = B.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues();
        Teuchos::ArrayRCP<Scalar>       Bdata = B.getValuesNonConst();
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
        TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::GESUM(C,alpha,A,beta,B,gamma): A and B must have the same dimensions.");
        Teuchos::RCP<Node> node = B.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues(),
                                        Bdata = B.getValues();
        Teuchos::ArrayRCP<Scalar>       Cdata = C.getValuesNonConst();
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

      inline static void Norm1(const MultiVector<Scalar,Node> &A, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        TEST_FOR_EXCEPTION(nC > Teuchos::as<size_t>(norms.size()), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Norm1(A,norms): norms must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          std::fill( norms.begin(), norms.begin() + nC, Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero() );
          return;
        }
        Teuchos::RCP<Node> node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues();
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
        Teuchos::RCP<Node> node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues(0);
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        SumAbsOp<Scalar> op;
        op.x = Adata(0,nR).getRawPtr();
        return node->parallel_reduce(0,nR,op);
      }

      inline static void Sum(const MultiVector<Scalar,Node> &A, const Teuchos::ArrayView<Scalar> &sums) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        TEST_FOR_EXCEPTION(nC > (size_t)sums.size(), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Sum(A,sums): sums must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          std::fill( sums.begin(), sums.begin() + nC, Teuchos::ScalarTraits<Scalar>::zero() );
          return;
        }
        Teuchos::RCP<Node> node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues();
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
        Teuchos::RCP<Node> node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues(0);
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
        Teuchos::RCP<Node> node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues(0);
        // prepare buffers
        ReadyBufferHelper<Node> rbh(node);
        rbh.begin();
        rbh.template addConstBuffer<Scalar>(Adata);
        rbh.end();
        MaxAbsOp<Scalar> op;
        op.x = Adata(0,nR).getRawPtr();
        return node->parallel_reduce(0,nR,op);
      }

      inline static void NormInf(const MultiVector<Scalar,Node> &A, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        TEST_FOR_EXCEPTION(nC > Teuchos::as<size_t>(norms.size()), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::NormInf(A,norms): norms must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          std::fill( norms.begin(), norms.begin() + nC, Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero() );
          return;
        }
        Teuchos::RCP<Node> node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues();
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

      inline static void Norm2Squared(const MultiVector<Scalar,Node> &A, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        TEST_FOR_EXCEPTION(nC > Teuchos::as<size_t>(norms.size()), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Norm2Squared(A,norms): norms must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          std::fill( norms.begin(), norms.begin() + nC, Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero() );
          return;
        }
        Teuchos::RCP<Node> node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues();
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
        Teuchos::RCP<Node> node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues(0);
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
        Teuchos::RCP<Node> node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues(0),
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

      inline static void WeightedNorm(const MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &weightVector, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride(),
                     Wstride = weightVector.getStride();
        TEST_FOR_EXCEPTION(nC > Teuchos::as<size_t>(norms.size()), std::runtime_error, 
            "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Norm1(A,norms): norms must have length as large as number of columns of A.");
        if (nR*nC == 0) {
          std::fill( norms.begin(), norms.begin() + nC, Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero() );
          return;
        }
        Teuchos::RCP<Node> node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues(),
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
        Teuchos::RCP<Node> node = A.getNode();
        Teuchos::ArrayRCP<Scalar> Adata = A.getValuesNonConst();
        // we'll overwrite all data covered by the multivector, but not off-stride data
        // therefore, we are write-only only in the case that stride=nR
        ReadWriteOption rw = (stride == nR ? WriteOnly : ReadWrite);
        Teuchos::ArrayRCP<Scalar> mvdata = node->template viewBufferNonConst<Scalar>(rw,stride*(nC-1)+nR,Adata);
        for (size_t j=0; j<nC; ++j) {
          for (size_t i=0; i<nR; ++i) {
            mvdata[j*stride + i] = SCT::random();
          }
        }
        mvdata = Teuchos::null;
      }

      inline static void Abs(MultiVector<Scalar,Node> &A, const MultiVector<Scalar,Node> &B) {
        const size_t nR = A.getNumRows();
        const size_t nC = A.getNumCols();
        const size_t Astride = A.getStride();
        const size_t Bstride = B.getStride();
        TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Abs(A,B): A and B must have the same dimensions.");
        if (nC*nR == 0) return;
        Teuchos::RCP<Node> node = A.getNode();
        Teuchos::ArrayRCP<const Scalar> Bdata = B.getValues();
        Teuchos::ArrayRCP<Scalar>       Adata = A.getValuesNonConst();
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
        TEST_FOR_EXCEPTION(nC != B.getNumCols() || nR != B.getNumRows(), std::runtime_error,
                           "DefaultArithmetic<" << Teuchos::typeName(A) << ">::Scale(B,alpha,A): A and B must have the same dimensions.");
        Teuchos::RCP<Node> node = B.getNode();
        Teuchos::ArrayRCP<const Scalar> Adata = A.getValues();
        Teuchos::ArrayRCP<Scalar>       Bdata = B.getValuesNonConst();
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
        Teuchos::RCP<Node> node = A.getNode();
        Teuchos::ArrayRCP<Scalar> data = A.getValuesNonConst();
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
                                   const Teuchos::ArrayRCP<Scalar> &values,
                                   size_t stride) {
        A.initializeValues(numRows,numCols,values,stride);
      }

      inline static Teuchos::ArrayRCP<const Scalar> getValues(const MultiVector<Scalar,Node> &A) {
        return A.getValues();
      }

      inline static Teuchos::ArrayRCP<const Scalar> getValues(const MultiVector<Scalar,Node> &A, size_t j) {
        return A.getValues(j);
      }

      inline static Teuchos::ArrayRCP<Scalar> getValuesNonConst(MultiVector<Scalar,Node> &A) {
        return A.getValuesNonConst();
      }

      inline static Teuchos::ArrayRCP<Scalar> getValuesNonConst(MultiVector<Scalar,Node> &A, size_t j) {
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

      inline static Teuchos::RCP<Node> getNode(const MultiVector<Scalar,Node> &A) {
        return A.getNode();
      }

  };

}

#endif
